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


#pragma once


#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

template <int dim>
class gp_t {

	static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6

	public:

	double strain_old[nvoi] = { 0.0 };
	double strain[nvoi];
	double stress[nvoi];
	double ctan[nvoi * nvoi];

	bool allocated; // flag for memory optimization

	double *vars_n; // vectors for calculations
	double *vars_k;
	double *u_n;
	double *u_k;
	int nvars;
	int nndim;

	long int cost;
	bool converged;
	bool subiterated;
	int coupling;

	gp_t():
		u_n(nullptr),
		u_k(nullptr),
		vars_n(nullptr),
		vars_k(nullptr),
		allocated(false),
		cost(0),
		converged(true),
		subiterated(false)
	{}

	~gp_t()
	{
		if (u_n != nullptr) {
			free(u_n);
		}
		if (u_k != nullptr) {
			free(u_k);
		}
		if (allocated) {
			free(vars_n);
			free(vars_k);
		}
	}

	void allocate_u()
	{
		u_n = (double *) calloc(nndim, sizeof(double));
		u_k = (double *) calloc(nndim, sizeof(double));
	}

	void allocate()
	{
		assert(!allocated);

		vars_n = (double *) calloc(nvars, sizeof(double));
		vars_k = (double *) calloc(nvars, sizeof(double));

		allocated = (vars_n && vars_k);
		assert(allocated);
	}


	void update_vars()
	{
		double *tmp = vars_n;
		vars_n = vars_k;
		vars_k = tmp;

		tmp = u_n;
		u_n = u_k;
		u_k = tmp;

		memcpy(strain_old, strain, nvoi * sizeof(double));
	}

	void write_restart(std::ofstream& file)
	{
		file.write((char *)&allocated, sizeof(bool));
		if (allocated) {
			file.write((char *)vars_n, nvars * sizeof(double));
			file.write((char *)u_n, nndim * sizeof(double));
		}
	}

	void read_restart(std::ifstream& file)
	{
		file.read((char *)&allocated, sizeof(bool));
		if (allocated) {
			vars_n = (double *)calloc(nvars, sizeof(double));
			vars_k = (double *)calloc(nvars, sizeof(double));

			file.read((char *)vars_n, nvars * sizeof(double));
			file.read((char *)u_n, nndim * sizeof(double));
		}
	}
};
