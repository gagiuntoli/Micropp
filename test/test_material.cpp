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
	// Execution ./test3d_1 n [print] [steps]
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [print = 0|1] [steps]" << endl;
		return(1);
	}

	const int dir = 1;
	const int n = atoi(argv[1]);

	const int print = (argc > 2) ? atoi(argv[2]) : 0;
	if (print < 0 || print > 1) {
		cerr << "Error in [print] argument, it only can be 0 or 1" << endl;
		return(1);
	}

	const int time_steps = (argc > 3 ? atoi(argv[3]) : 10);
	const int micro_type = MIC_SPHERES;
	const double micro_params[4] = { 1.0, 1.0, 1.0, 0.2 };
	const int size[3] = { n, n, n };

	ofstream file;
	file.open("result.dat");

	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.3, 5.0e4, 2.0e4, 1);
	mat_params[1].set(1.0e3, 0.3, 5.0e4, 1.0e3, 1);

	material_damage mat_damage(1.0e6, 0.3, 5.0e4);

	//micropp<3> micro(1, size, micro_type, micro_params, mat_params, NO_COUPLING);
	micropp<3> micro(1, size, micro_type, micro_params, mat_params, ONE_WAY, true, 5);
	micro.print_info();

	file.close();
	return 0;
}
