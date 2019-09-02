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

#include "micropp.hpp"

using namespace std;

template <int tdim>
class test_t : public micropp<tdim> {
	public:
		test_t(const micropp_params_t params)
			:micropp<tdim> (params)
		{};

		~test_t() {};

		void public_get_elem_nodes(int n[8], int ex, int ey, int ez = 0)
		{
			micropp<tdim>::get_elem_nodes(n, ex, ey, ez);
		};
};

int main (int argc, char *argv[])
{
	int n[8];

	micropp_params_t mic_params;
	mic_params.size[0] = 5;
	mic_params.size[1] = 5;
	mic_params.size[2] = 5;

	test_t<3> test(mic_params);

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
