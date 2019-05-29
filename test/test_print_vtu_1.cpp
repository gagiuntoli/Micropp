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
		cerr << "Usage: " << argv[0] << " <n> <mic_type = [0 ... 8]>" << endl;
		exit(1);
	}

	const int n = atoi(argv[1]);
	const int size[3] = { n, n, n };

	const int mic_selected = atoi(argv[2]);
	if (mic_selected < 0 || mic_selected > MIC3D_8) {
		cerr << "<mic_type = [0 ... 8]>" << endl;
		exit(1);
	}

	material_base materials[3];
	material_set(&materials[0], 0, 1.0e6, 0.3, 5.0e4, 2.0e4, 0.0);
	material_set(&materials[1], 1, 1.0e3, 0.3, 5.0e4, 1.0e3, 0.0);
	material_set(&materials[2], 1, 1.0e3, 0.3, 5.0e4, 1.0e3, 0.0);

	double params[9][4] = {
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.1 },
		{ 1.0, 1.0, 1.0, 0.6 },
		{ 1.0, 1.0, 1.0, 0.1 }
	};

	micropp<3> *micro = new micropp<3>(1, size, mic_selected, params[mic_selected],
					   materials, nullptr, false, 0, 0,
					   NR_MAX_ITS, NR_MAX_TOL, NR_REL_TOL,
					   false);
	micro->print_info();

	char filename[128];
	snprintf(filename, 128, "micro_type_%d", mic_selected);
	micro->output (0, filename);
	delete micro;

	return 0;
}
