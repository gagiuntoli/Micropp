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

#ifndef UTIL_HPP
#define UTIL_HPP

// Debug print macro.
#ifdef NDEBUG
#define dprintf(...)
#else
#define dprintf(...) fprintf(stderr, __VA_ARGS__)
#endif


#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>


#include <cmath>
#include <ctime>
#include <cstring>
#include <cassert>


using namespace std;


inline uint64_t devest(const vector<uint64_t>  &in, const uint64_t mean)
{
	uint64_t out = 0;
	for (const auto &x : in) {
		const uint64_t tmp = (x - mean);
		out += tmp * tmp;
	}

	return sqrt(out / in.size());
}


inline void filter(double *arr, int n, double rel_tol)
{
#ifdef FILTER
	double max = arr[0];
	for (int i = 1; i < n; ++i)
		if (arr[i] > max)
			max = arr[i];

	for (int i = 0; i < n; ++i)
		arr[i] = (fabs(arr[i]) > max * rel_tol) ? arr[i] : 0.0;
#endif
}


constexpr int mypow(int v, int e)
{
	return (e == 0) ? 1 : v * mypow(v, e - 1);
}


inline void print_vec(const double *vec, int n, const char file_name[])
{
	FILE *file = fopen(file_name, "w");
	for (int i = 0; i < n; ++i)
		fprintf(file, "[%lf]\n", vec[i]);
	fclose(file);
}


inline void mvp_2(const double m[2][2], const double x[2], double y[2])
{
	for (int i = 0; i < 2; ++i) {
		double tmp = 0.0;
		for (int j = 0; j < 2; ++j)
			tmp += m[i][j] * x[j];
		y[i] = tmp;
	}
}


inline void mvp_3(const double m[2][2], const double x[2], double y[2])
{
	for (int i = 0; i < 2; ++i) {
		double tmp = 0.0;
		for (int j = 0; j < 2; ++j)
			tmp += m[i][j] * x[j];
		y[i] = tmp;
	}
}


inline void mvp_3(const double m[3][3], const double x[3], double y[3])
{
	for (int i = 0; i < 3; ++i) {
		double tmp = 0.0;
		for (int j = 0; j < 3; ++j)
			tmp += m[i][j] * x[j];
		y[i] = tmp;
	}
}


#endif
