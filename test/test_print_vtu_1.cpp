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
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " n [mic_type = [0 ... 7]]" << endl;
		return(1);
	}

	const int n = atoi(argv[1]);
	const int size[3] = { n, n, n };

	const int mic_selected = (argc == 3) ? atoi(argv[2]) : -1;
	if (mic_selected > MIC_SPHERES) {
		cerr << "mic_type > " << MIC_SPHERES << endl;
		return(1);
	}


	material_base materials[2];
	material_set(&materials[0], 0, 1.0e6, 0.3, 5.0e4, 2.0e4, 0.0);
	material_set(&materials[1], 1, 1.0e3, 0.3, 5.0e4, 1.0e3, 0.0);

	double params[8][4] = {
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.6 }
	};

	if (mic_selected < 0) {

		for(int micro_type = 0; micro_type < MIC_SPHERES + 1; ++micro_type) {

			cout << "Plotting Micro Type : " << micro_type << endl;
			micropp<3> *micro = new micropp<3>(1, size, micro_type,
							   params[micro_type],
							   materials, NO_COUPLING);
			micro->print_info();

			char filename[128];
			snprintf(filename, 128, "micro_type_%d", micro_type);
			micro->output (0, filename);
			delete micro;
		}

	} else {

		cout << "Plotting Micro Type : " << mic_selected << endl;
		micropp<3> *micro = new micropp<3>(1, size, mic_selected,
						   params[mic_selected],
						   materials, NO_COUPLING);
		micro->print_info();

		char filename[128];
		snprintf(filename, 128, "micro_type_%d", mic_selected);
		micro->output (0, filename);
		delete micro;
	}


	return 0;
}
