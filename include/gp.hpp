/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                        Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include <cassert>
#include <cstdlib>

#include "newton.hpp"

template <int dim>
class gp_t {

	static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6

	public:

	double macro_strain[nvoi];
	double macro_stress[nvoi];
	double *macro_ctan;

	bool allocated; // flag for memory optimization

	double *int_vars_n; // vectors for calculations
	double *int_vars_k;
	double *u_n;
	double *u_k;

	newton_t newton;
	long int cost;

	gp_t():
		u_n(nullptr),
		u_k(nullptr),
		int_vars_n(nullptr),
		int_vars_k(nullptr),
		allocated(false)
	{}

	~gp_t()
	{
		free(u_n);
		free(u_k);
		if (allocated) {
			free(macro_ctan);
			free(int_vars_n);
			free(int_vars_k);
		}
	}

	void allocate(const int num_int_vars)
	{
		assert(!allocated);

		macro_ctan = (double *) calloc(nvoi * nvoi, sizeof(double));
		int_vars_n = (double *) calloc(num_int_vars, sizeof(double));
		int_vars_k = (double *) calloc(num_int_vars, sizeof(double));

		allocated = (int_vars_n && int_vars_k);
		assert(allocated);
	}


	void update_vars()
	{
		double *tmp = int_vars_n;
		int_vars_n = int_vars_k;
		int_vars_k = tmp;

		tmp = u_n;
		u_n = u_k;
		u_k = tmp;
	}
};
