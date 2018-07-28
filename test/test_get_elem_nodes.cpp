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
		test_t(const int size[3], const double micro_params[5],
		       const material_t mat_params[2])
			:micropp<tdim> (4, size, 1, micro_params, mat_params)
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
	const double micro_params[5] = { 1., 1., 1., 0.1, 0. };

	int n[8];

	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.3, 5.0e4, 5.0e4, 1);
	mat_params[1].set(1.0e6, 0.3, 1.0e4, 0.0e-1, 0);

	test_t<3> test(size, micro_params, mat_params);

	test.public_get_elem_nodes(n, 0, 0, 0);
	const int n_1_exact[8] = { 0, 1, 6, 5, 25, 26, 31, 30 };

	for (int i = 0; i < 8; ++i) {
		printf("Asserting: n[%d] = %d == %d\n", i, n[i], n_1_exact[i]);
		assert(n[i] == n_1_exact[i]);
	}

	test.public_get_elem_nodes(n, 0, 1, 0);
	const int n_2_exact[8] = { 5, 6, 11, 10, 30, 31, 36, 35 };

	for (int i = 0; i < 8; ++i) {
		printf("Asserting: n[%d] = %d == %d\n", i, n[i], n_2_exact[i]);
		assert(n[i] == n_2_exact[i]);
	}

	test.public_get_elem_nodes(n, 0, 1);
	const int n_3_exact[8] = { 5, 6, 11, 10, 30, 31, 36, 35 };

	for (int i = 0; i < 8; ++i) {
		printf("Asserting: n[%d] = %d == %d\n", i, n[i], n_3_exact[i]);
		assert(n[i] == n_3_exact[i]);
	}

	return 0;
}
