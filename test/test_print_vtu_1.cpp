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
		cerr << "Usage: " << argv[0] << " nx ny nz" << endl;
		return(1);
	}
	const int size[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };

	const double special_param[5] = { .11, .11, .11, .11, .15 };
	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.3, 5.0e4, 5.0e4, 1);
	mat_params[1].set(1.0e6, 0.3, 1.0e4, 0.0e-1, 0);

	for(int micro_type = 0; micro_type < 6; ++micro_type) {

		cout << "Plotting Micro Type : " << micro_type << endl;
		double micro_params[4] = { 1.0, 1.0, 1.0, special_param[micro_type] };
		micropp<3> *micro = new micropp<3>(1, size, micro_type,
						   micro_params, mat_params);


		char filename[128];
		snprintf(filename, 128, "micro_type_%d", micro_type);
		micro->output (0, filename);
		delete micro;
	}

	return 0;
}
