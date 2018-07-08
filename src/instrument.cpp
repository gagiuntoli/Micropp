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

#ifndef NDEBUG // For optimized compilation

#include "instrument.hpp"
#include <iomanip>      // std::setw

using namespace std;

unordered_map<string, timevect> instrument::times;                   // Final register time
thread_local unordered_map<string, timestack> instrument::startTime; // function stack calls
struct timespec instrument::initialTime = {0};
atomic<int> instrument::instances(0);

void instrument::init()
{
	if (instances++)
		getTime(initialTime);
}

void instrument::finalize()
{
	if (0 == --instances) {
		struct timespec finalTime;
		getTime(initialTime);

		size_t cont = 0;
		cout << "# Final execution report" << endl;

		for (auto const &f : times) {
			const double total = accumulate(f.second.begin(), f.second.end(),
			                                0.0,
			                                [](double x, const struct timespec &ts) {
				                                return x + getus(ts);
			                                });
			const size_t entries = f.second.size();
			const double mean = total / entries;

			cout << cont++ << " "
			     << setw(20) << left
			     << f.first << " "
			     << entries << " "
			     << total << " "
			     << mean << " "
			     << devest(f.second, mean)
			     << endl;

			if (!startTime[f.first].empty())
				cout << "# " << startTime[f.first].size()
				     << " Unpaired entries for: " << f.first << endl;
		}

		cout << "Total Time: " << getus(getDt(initialTime, finalTime)) << endl;
	}
}

#endif // NDEBUG
