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

template <int tdim>
class test_t : public micropp<tdim> {

	public:
		test_t(const int size[3], const double mat_params[20],
		       const int micro_types[2], const double params[5])
			: micropp<tdim> (4, size, 1, mat_params, micro_types, params)
		{};

		~test_t() {};

        void public_get_elem_nodes(int n[8], int ex, int ey, int ez = 0)
		{
			micropp<tdim>::get_elem_nodes(n, ex, ey, ez);
		};
};

int main (int argc, char *argv[])
{
	const int size[3] = { 5, 5, 5 };
	double mat_params[20];
	mat_params[0 * MAX_MAT_PARAM + 0] = 1.0e6; // E
	mat_params[0 * MAX_MAT_PARAM + 1] = 0.3;   // nu
	mat_params[0 * MAX_MAT_PARAM + 2] = 5.0e4; // Sy
	mat_params[0 * MAX_MAT_PARAM + 3] = 5.0e4; // Ka

	mat_params[1 * MAX_MAT_PARAM + 0] = 1.0e6;
	mat_params[1 * MAX_MAT_PARAM + 1] = 0.3;
	mat_params[1 * MAX_MAT_PARAM + 2] = 1.0e4;
	mat_params[1 * MAX_MAT_PARAM + 3] = 0.0e-1;

	const int micro_types[2] = { 1, 0 };
	const double params[5] = { 1., 1., 1., .1, 0. };

	test_t<2> test(size, mat_params, micro_types, params);

	int n[8];
	test.public_get_elem_nodes(n, 0, 0, 0);

	return 0;
}
