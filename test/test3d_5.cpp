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
		cerr << "Usage: " << argv[0] << " nx ny nz [steps]" << endl;
		return(1);
	}

	const int nx = atoi(argv[1]);
	const int ny = atoi(argv[2]);
	const int nz = atoi(argv[3]);
	const int time_steps = (argc > 4 ? atoi(argv[4]) : 10);  // Optional value
	const int ngp = 10;

	assert(nx > 1 && ny > 1 && nz > 1 && time_steps > 0);

	int size[dim] = {nx, ny, nz};

	int micro_type = 1;	// 2 materiales matriz y fibra (3D esfera en matriz)

	double micro_params[5] = {1.0,		// lx
		1.0,		// ly
		1.0,		// lz
		0.1,		// layer width
		1.0e-5};	// inv_max

	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.3, 5.0e4, 5.0e4, 1);
	mat_params[1].set(1.0e6, 0.3, 1.0e4, 0.0e-1, 0);

	int dir = 2;
	const double d_eps_1 = 0.01, d_eps_2 = -0.008;;
	double eps_1[6] = {0.0}, eps_2[6] = {0.0};
	double sig[6], (*sig_test)[3];

	sig_test = (double (*)[3]) malloc (3 * ngp * sizeof(double));

	micropp<dim> micro(ngp, size, micro_type, micro_params, mat_params);

	for (int t = 0; t < time_steps; ++t) {
		cout << "Time step = " << t << endl;
		if (t < 30) {
			eps_1[dir] += d_eps_1;
			eps_2[dir] += d_eps_2;
		} else if (t < 80) {
			eps_1[dir] -= d_eps_1;
			eps_2[dir] -= d_eps_2;
		} else if (t < 130) {
			eps_1[dir] += d_eps_1;
			eps_2[dir] += d_eps_2;
		}

		cout << "setting strains ..." << endl;
		for (int gp = 0; gp < ngp; ++gp) {
			micro.set_macro_strain(gp, (gp%2) ? eps_1 : eps_2);
			cout << "gp = " << gp << " eps = ";
			cout << scientific;

			for (int i = 0; i < 6; ++i)
				cout << setw(14) << ((gp%2) ? eps_2[i] : eps_1[i])<< " ";
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

			memcpy(sig_test[gp], sig, 3 * sizeof(double));
		}

		for (int gp = 2; gp < ngp; ++gp) {
			cout << "Diff: \t";
			for (int i = 0; i < 3; ++i) {
				const double tmp = fabs(sig_test[gp][i] - sig_test[gp%2 ? 1 : 0][i]);
				cout << tmp << "\t";
				assert(tmp < 1.0e-6);
			}
			cout << endl;
		}

		micro.update_vars();

		char filename[128];
		snprintf(filename, 128, "micropp_%d", t); 
		micro.output (1, filename);
		cout << endl;
	}
	return 0;
}
