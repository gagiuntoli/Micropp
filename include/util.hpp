/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
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

inline double getus(const struct timespec &ts)
{
	return ts.tv_sec * 1000000.0 + double(ts.tv_nsec) * 1.0E-3;
}

inline double devest(const vector<struct timespec>  &in, double mean)
{
	double out = 0;
	for (const auto &x : in) {
		const double tmp = (getus(x) - mean);
		out += tmp * tmp;
	}

	out /= in.size();

	return sqrt(out);
}

inline void getTime(struct timespec &ts)
{
	const int rc = clock_gettime(CLOCK_MONOTONIC, &ts);
	if (rc) {
		std::cerr << "Error reading time: " << strerror(errno) << std::endl;
		exit(1);
	}
}

static inline struct timespec getDt(const struct timespec &t1,
                                    const struct timespec &t2)
{
	struct timespec dt;

	if (t2.tv_nsec >= t1.tv_nsec) {
		dt.tv_sec = t2.tv_sec - t1.tv_sec;
		dt.tv_nsec = t2.tv_nsec - t1.tv_nsec;
	} else {
		dt.tv_sec = t2.tv_sec - 1 - t1.tv_sec;
		dt.tv_nsec = 1000000000L + t2.tv_nsec - t1.tv_nsec;
	}
	return dt;
}

#endif //UTIL_HPP
