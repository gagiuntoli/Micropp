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

#define D_EPS 5.0e-4

int main (int argc, char *argv[])
{
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [steps]" << endl;
		return(1);
	}

	const int dir = 1;
	const int n = atoi(argv[1]);
	const int time_steps = (argc > 2 ? atoi(argv[2]) : 10);  // Optional value
	const int micro_type = MIC_LAYER_Y;

	ofstream file;
	file.open("result.dat");

	micropp<3> *micropp3 = get_instance_3_simple (1, n, micro_type, 1.0);
	micropp3->print_info();

	double sig[6];
	double eps[6] = { 0. };

	for (int t = 0; t < time_steps; ++t) {

		cout << "time step = " << t << endl;

		if (t < 80)
			eps[dir] += D_EPS;
		else if (t < 160)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		micropp3->set_macro_strain(0, eps);
		micropp3->homogenize();
		micropp3->get_macro_stress(0, sig);
		int newton_its = micropp3->get_sigma_newton_its(0);
		int non_linear = micropp3->is_non_linear(0);

		char filename[128];
		snprintf(filename, 128, "micro_type_%d", micro_type);
		micropp3->output (0, filename);

		micropp3->update_vars();

		cout << "non_linear = \t" << non_linear << "\tnewton its =\t" << newton_its << endl;
		cout << "eps =\t";
		for (int i = 0; i < 6; ++i)
			cout << setw(14) << eps[i] << "\t";
		cout << endl;

		cout << "sig =\t";
		for (int i = 0; i < 6; ++i)
			cout << setw(14) << sig[i] << "\t";
		cout << endl;

		cout << endl;
		file << setw(14) << eps[dir] << "\t" << sig[dir] << "\t" << endl;
	}

	file.close();
	return 0;
}
