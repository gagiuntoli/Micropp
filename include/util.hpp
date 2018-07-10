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

#endif //UTIL_HPP
