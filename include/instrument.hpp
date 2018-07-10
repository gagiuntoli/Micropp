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


#ifndef INSTRUMENT_HPP
#define INSTRUMENT_HPP

#ifdef NDEBUG // For optimized compilation

#define INST_CONSTRUCT
#define INST_DESTRUCT
#define INST_START

#else // For debug compilation

#define INST_CONSTRUCT instrument::initialize()
#define INST_DESTRUCT instrument::finalize()
#define INST_START instrument __timer__(__FUNCTION__)

#include <chrono>
#include <unordered_map>
#include <vector>
#include <stack>
#include <string>
#include <atomic>

#include "util.hpp"

using namespace std;

typedef vector<uint64_t> timevect;

class instrument {
	private:
		const uint64_t start_time_;
		const string funct_;

		static atomic<size_t> instances;                     // Counter
		static uint64_t initialTime;                         // Time collection
		static unordered_map<string, timevect> times;        // Final register time

		static inline uint64_t take_time_stamp()
		{
			return uint64_t(
				chrono::high_resolution_clock::now().time_since_epoch().count());
		}

	public:

		instrument(const string funct):
			funct_(funct), start_time_(take_time_stamp())
		{}

		~instrument()
		{
			const uint64_t elapsed =  (take_time_stamp() - start_time_) * 1E-3;
			times[funct_].push_back(elapsed);
		}

		static void initialize();

		static void finalize();
};

#endif // NDEBUG

#endif //INSTRUMENT_HPP

