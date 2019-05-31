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

	const int ngp = 10;
	const int time_steps = (argc > 4 ? atoi(argv[4]) : 10);  // Optional value

	micropp_params_t mic_params;

	mic_params.ngp = ngp;
	mic_params.size[0] = atoi(argv[1]);
	mic_params.size[1] = atoi(argv[2]);
	mic_params.size[2] = atoi(argv[3]);
	mic_params.type = MIC_SPHERE;
	material_set(&mic_params.materials[0], 0, 1.0e6, 0.3, 5.0e4, 2.0e4, 0.0);
	material_set(&mic_params.materials[1], 1, 1.0e3, 0.3, 5.0e4, 1.0e3, 0.0);
	material_set(&mic_params.materials[2], 0, 1.0e3, 0.3, 5.0e4, 1.0e3, 0.0);

	micropp<dim> micro(mic_params);

	int dir = 2;
	const double d_eps_1 = 0.01, d_eps_2 = -0.008;;
	double eps_1[6] = {0.0}, eps_2[6] = {0.0};
	double sig[6], (*sig_test)[3];

	sig_test = (double (*)[3]) malloc (3 * ngp * sizeof(double));

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
			micro.set_strain(gp, (gp%2) ? eps_1 : eps_2);
			cout << "gp = " << gp << " eps = ";
			cout << scientific;

			for (int i = 0; i < 6; ++i)
				cout << setw(14) << ((gp%2) ? eps_2[i] : eps_1[i])<< " ";
			cout << endl;
		}

		cout << "Homogenizing ..." << endl;
		micro.homogenize();

		cout << "Getting stresses ..." << endl;
		for (int gp = 0; gp < ngp; ++gp) {
			micro.get_stress(gp, sig);
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
