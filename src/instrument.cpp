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

#ifdef TIMER // For optimized compilation

#include "instrument.hpp"
#include <iomanip>      // std::setw

using namespace std;

atomic<size_t> instrument::instances(0);
uint64_t instrument::initialTime(0);
unordered_map<string, timevect> instrument::times; // Final register time

void instrument::initialize()
{
	if (!instances++)
		initialTime = take_time_stamp();
}

void instrument::finalize()
{
	if (0 == --instances) {
		const uint64_t elapsed = (take_time_stamp() - initialTime) * 1E-3;

		size_t cont = 0;
		cout << "# Final execution report: total time = " << elapsed
			<< endl;
		cout << setw(6) << left << "#No"
			<< setw(25) << "function"
			<< setw(8) << right << "calls"
			<< setw(16) << "total time"
			<< setw(16) << "percent"
			<< setw(16) << "mean"
			<< setw(16) << "stdev"
			<< setw(16) << "relative"
			<< endl;

		cout.precision(2);
		cout << fixed;

		for (auto const &f : times) {
			const uint64_t total = accumulate(f.second.begin(),
							  f.second.end(),
							  0.0);
			const size_t entries = f.second.size();
			const uint64_t mean = total / entries;
			const double percent = double(total) * 100 / elapsed;
			const uint64_t stdev = devest(f.second, mean);
			const double relative = double(stdev) * 100 / mean;

			cout << setw(6) << left << cont++
				<< setw(25) << f.first
				<< setw(8) << right << entries
				<< setw(16) << total
				<< setw(16) << percent
				<< setw(16) << mean
				<< setw(16) << stdev
				<< setw(16) << relative
				<< endl;
		}
	}

}

#endif // TIMER
